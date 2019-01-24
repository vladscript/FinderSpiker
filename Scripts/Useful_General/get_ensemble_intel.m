% Get Data from Several Analysis Variables
function [Rasters,LENGHTRASTER,NGroups,SIG_FRAMES,LABELS,THR,...
    ExperimentRaster,signif_frames,labels_frames,CummFrames,CummGroups,...
    TotalNG]=get_ensemble_intel(NN_data_Read)
NC=numel(NN_data_Read);
% DATA PER CONDITION
Rasters={};                                 % Analyzed Raster
LENGHTRASTER={};                            % Frames per Condition
NGroups={};                                 % N ensembles per Condition
SIG_FRAMES={};                              % Frames above CAG Threshold
LABELS={};                                  % Ensemble Labels of the Frames
THR={};                                     % CAG threshold per condition
% CONCATENATED DATA
ExperimentRaster=[];    % Matrix Rx to concatenate Rasters
signif_frames=[];       % Frames with Significant Coactivity
labels_frames=[];       % Ensemble Labels
CummFrames=0;           % Cummulative Frames
CummGroups=0;           % Cummulative Groups
TotalNG=0;              % Total Nuber of Groups (N conditions->NG ensembles)

% Main Loop
for c=1:NC
    NN_data=NN_data_Read{c};
    % Get RASTER data from Structure
    Rasters{c}=NN_data.Data.Data;                       % Raster: Frames x Cells
    LENGHTRASTER{c}=length(NN_data.Data.Data);          % Raster's Lengths
    NamesFields=fieldnames(NN_data);
    % Get ENSEMBLES data from Structure
    if ismember('Clustering',NamesFields)   % Ensembles Activity
        NGroups{c}=NN_data.Clustering.TotalStates;          % Number of Groups
        SIG_FRAMES{c}=find(NN_data.Peaks.Index);            % Frames with CAG above Threshold
        LABELS{c}=NN_data.Clustering.VectorStateIndex+CummGroups;      % Ensemble Labels
        THR{c,1}=NN_data.Peaks.Threshold;                   % CAG Threshold
        fprintf('>> Ensembles in Condition %i \n',c)
        tabulate(LABELS{c});
    else                                    %  NO Ensembles Activity
        fprintf('>> No Ensembles in Condition %i \n',c)
        NGroups{c}=0;      % Number of Groups
        SIG_FRAMES{c}=[];  % Frames with CAG above Threshold
        LABELS{c}=[];      % Ensemble Labels
    end
    % Concatenation **************
    % Raster
    ExperimentRaster=[ExperimentRaster;Rasters{c}]; % Frames x Cells
    % Significative Frames
    signif_frames=[makerowvector(signif_frames), makerowvector(SIG_FRAMES{c})+CummFrames];
    % Total Number of Frames
    CummFrames=CummFrames+LENGHTRASTER{c};
    % Frame Labels
    labels_frames=[labels_frames;LABELS{c}];
    % Total Number of Groups
    CummGroups=CummGroups+NGroups{c};
    % Group COunter:
    TotalNG=TotalNG+NGroups{c};
end