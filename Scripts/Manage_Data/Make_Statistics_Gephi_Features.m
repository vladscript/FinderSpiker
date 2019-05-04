%% Make Data set of Statistics of Netowrk Parameters
%% 0. Load MAT File
[FileName,PathName] = uigetfile('.mat','Load DATASET of GEPHI Features');
fullFileName = fullfile(PathName, FileName);
load(fullFileName);
% setup
[NE,NC]=size(GEPHIDATA); % Number of Experiments & Conditions
%% 1. Chooose Feature from table: clustering, grade, etc, etc
NetFeaturesNames=GEPHIDATA{1}.Properties.VariableNames(6:end)';
[index_var,index_CHECK] = listdlg('PromptString','Select a Network Feature:',...
                    'SelectionMode','single',...
                    'ListString',NetFeaturesNames);
ActualFeature=NetFeaturesNames(index_var);
% 1.1 Check if the features is in ALL tables
ExsitFeature=zeros(NE,NC);
ColumnIndx=zeros(NE,NC);
for e=1:NE
    for c=NC
        GEPHIDATA{e,c}
        ExsitFeature(e,c)
    end
end

% 2. Check Experiment ID's are OK (non repeated and matched)
% 3. Read vector from ALL Tables and ALL Conditions
% 4. Make Statistics: mean, variance, mode, median, etc
% 5. Make Table for ML o Statistical Anlysis
