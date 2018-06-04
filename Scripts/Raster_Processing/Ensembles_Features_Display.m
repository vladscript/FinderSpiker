%% Raster Ensembles Display
% Script that read and plot Ensembles Features
% between or among Conditions of one Several Experiments
% 'Active Neurons','Duration','MeanActivity','EffectiveActivity',...
% 'ISImean','ISImode','ISIvar','ISIskew','ISIkurt',...
% 'Lengthmean','Lengthmode','Lengthvar','Lengthskew','Lengthkurt'
%% Read CSV Files
NC = inputdlg('Number of Conditions: ',...
             'Raster Features', [1 75]);
NC = str2num(NC{:});    
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
        % Only Dominant ensmeble Intel
        ENSDOM=rowFeatures(5);
        Nens=rowFeatures(1);
        rowFeaturesOK=rowFeatures(1:8);
        rowFeaturesOK(9)=rowFeatures(8+ENSDOM);
        rowFeaturesOK(10)=rowFeatures(8+Nens+ENSDOM);
        FeaturesSingle=[FeaturesSingle;rowFeaturesOK];
        
    end
    ENSEMBLE_FEATURES{i}=FeaturesSingle;
    ENSEMBLE_NAMES{i}=Ensemble_Names;
end
%% Plot Data
% First Figure:
%   N-Ensembles,Dunns Index,RateOfTransitions,RateOfCycles
FeaturesA=figure;
FeaturesA.Name='Ensemble Features';
h1=subplot(2,2,1); % N-Ensembles
plot_box(ENSEMBLE_NAMES,ENSEMBLE_FEATURES,Names_Conditions,1,h1)
title(h1,'N-ensembles')
h2=subplot(2,2,2); % Dunn's Index
plot_box(ENSEMBLE_NAMES,ENSEMBLE_FEATURES,Names_Conditions,2,h2)
title(h2,'Index Dunnn')
h3=subplot(2,2,3); % Rate of Transitions [min]
plot_box(ENSEMBLE_NAMES,ENSEMBLE_FEATURES,Names_Conditions,3,h3)
title(h3,'Rate of Transition [T/min]')
h4=subplot(2,2,4); % Rate of Cycles per MINUTE
plot_box(ENSEMBLE_NAMES,ENSEMBLE_FEATURES,Names_Conditions,4,h4)
title(h4,'Rate of Cycles [C/min]')

%% Pie Charts of Kind of Cycles
% 'ISImean','ISImode','ISIvar','ISIskew','ISIkurt',...
% 'Lengthmean','Lengthmode','Lengthvar','Lengthskew','Lengthkurt'
FeaturesB=figure;
FeaturesB.Name='Type of Cylces or Paths';
% [Simple|Open|Closed]
ColorPie=[44,188,187;110,59,128;165,71,11]/255;
for c=1:NC
    ax=subplot(1,NC,c);
    F=ENSEMBLE_FEATURES{c}(:,6:8);
    [Nr,~]=size(F);
    for r=1:Nr
        F(r,isnan(F(r,:)))=0;
        MINF=min(F(r,F(r,:)>0));
        if isempty(MINF)
            MINF=0;
        end
        if MINF==0
            N=0;
        else
            N=1/MINF;
        end
        auxF(r,:)=N*F(r,:);
    end
    Pf=pie(sum(auxF));
    Pf(1).FaceColor=ColorPie(1,:);
    Pf(3).FaceColor=ColorPie(2,:);
    Pf(5).FaceColor=ColorPie(3,:);
    title(ax,Names_Conditions{c})
    l=legend('Simple','Open','Closed');
    l.Location='southoutside';
    l.Box='off';
end

%% Dominant Ensemble Activity
FeaturesC=figure;
FeaturesC.Name='Dominant Ensemble Features';
g1=subplot(1,2,1); % ISImode
plot_box(ENSEMBLE_NAMES,ENSEMBLE_FEATURES,Names_Conditions,9,g1)
title(g1,'Rate')

g2=subplot(1,2,2); % ISImode
plot_box(ENSEMBLE_NAMES,ENSEMBLE_FEATURES,Names_Conditions,10,g2)
title(g2,'Dominance Index')

