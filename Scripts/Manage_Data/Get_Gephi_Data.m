% Function to accumulate Tables Generated in Gephi
% of each Experiment for different Conditions

% Input
%   N Conditions: by GUI
%   Load experiment manually from different folder: one by one
%   
% Output varbales @ mat Files
%    File .mat with all the Network features per Experiment
%% Setup
% Initial:
clear; close all; clc;
runs=1;             % Runs Counter
EXPS={};            % List Of Experiments

% N Conditions:
NC = inputdlg('Number of Conditions: ',...
             'Network Gephi Features', [1 75]);
NC = str2double(NC{:});    

% Name Conditions
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

% Directory:
Dirpwd=pwd;
slashesindx=find(Dirpwd=='\');
CurrentPathOK=[Dirpwd(1:slashesindx(end))];
% Load File 
c=1; % naive initialization
[FileName,PathName,MoreFiles] = uigetfile('*.csv',['Condition: ',num2str(c),' Gephi Data .csv file'],...
    'MultiSelect', 'off',CurrentPathOK);
%% Collect DATA
GEPHIDATA={}; % Data set of Network Features of ALl Experiments
for c=1:NC
    fprintf('>>Loading Data from Ccndition %i of %i:\n',c,NC);
    auxc=1;
    while MoreFiles
        X=readtable([PathName,FileName]);
        [Nnodes,NFeatures]=size(X);
        SpaceIndx=find(FileName==' ',1);
        EXPID=FileName(1:SpaceIndx-1);
        Y=[table(repmat(Names_Conditions{c},Nnodes,1),'VariableNames',{'Condition'}),...
            table(repmat(EXPID,Nnodes,1),'VariableNames',{'EXP_ID'})...
            X];
        GEPHIDATA{auxc,c}=Y;
        % Disp Datasets Selected:
        EXPS{runs,1}=FileName
        CurrentPathOK=PathName;
        runs=runs+1;
        [FileName,PathName,MoreFiles] = uigetfile('*.csv',['Condition: ',num2str(c),...
            ' : ',Names_Conditions{c},' Gephi Data .csv file| Press CANCEL to start the Following Condition'],...
        'MultiSelect', 'off',CurrentPathOK);
        auxc=auxc+1;
    end
    % To get into the (fruit) loop again
    if c<NC
        MoreFiles=true;
        [FileName,PathName,MoreFiles] = uigetfile('*.csv',['Condition: ',num2str(c+1),...
            ' : ',Names_Conditions{c+1},' Gephi Data .csv file'],...
        'MultiSelect', 'off',CurrentPathOK);
    end
end
disp('>>Loading Gephi Data :Done.')
disp('>>Saving Gephi Data ...')
%% Save  Data

% Input dialogue
CurrentPathOK=[Dirpwd(1:slashesindx(end))];

[file,path] = uiputfile('myNetwrokDataSet.mat','Save file name');
if file ~=0
    fullFileName = fullfile(path, file);
    save(fullFileName,'GEPHIDATA','Names_Conditions');
    disp('>>Network Dataset: Saved.')
else
    disp('>>CANCEL')
end
%% FINISH
clear; close all; clc; disp('>>DONE')