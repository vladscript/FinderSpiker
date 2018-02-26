%% Function to Plot Ensembles **************************************
% Input
%   R_Condition         Rasters by Condition
%   EnsembleName        Condition Names
%   Ensembled_Labels    Labels in time
%   Ensemble_Threshold  Significane Threshold
%   ColorState          Color Map 
%   fs                  Sampling Frequency
% Output
%   Figure of Ensembles Whole Experiment
function Plot_Ensembles_Experiment(R_Condition,EnsembleName,...
    Ensembled_Labels,Ensemble_Threshold,UniRMutiE,ColorState,fs)
%% Set up: Concatenate Stuff
NC=numel(R_Condition);  % N Conditions
NE=0;                   % N Ensembles Counter
ExperimentRaster=[];    % Whole Raster
labels_frames=[];       % Ensemble Labels
signif_frames=[];       % Frames with Significant Coactivity
CummFrames=0;           % Cummulative Frames
for i=1:NC
    % Raster
    ExperimentRaster=[ExperimentRaster,R_Condition{i}]; % Cells x Frames
    % Significative Frames
    signif_frames=[signif_frames, find(sum(R_Condition{i})>=Ensemble_Threshold{i})+CummFrames];
    CummFrames = CummFrames + length( R_Condition{i} );
    NG=numel(unique(Ensembled_Labels{i}));           % N ensembles @ condition 
    % Frame Labels
    if UniRMutiE
        labels_frames=[labels_frames; Ensembled_Labels{i}];
    else
        labels_frames=[labels_frames; Ensembled_Labels{i}+NE];
    end
    NE=NE+NG;
    % Threshold
    THR=Ensemble_Threshold{i};
end

[Cells,~]=size(ExperimentRaster);
CoAc=sum(ExperimentRaster);
Indexes=1:Cells;
[New_Order_Clustering,~]=OrderClusters(labels_frames,signif_frames,ExperimentRaster',NE);
%% PLOT
%   Original *************************************************************
OriginalExperiment=ExperimentRaster';
Plot_Raster_Ensembles(OriginalExperiment,Indexes,5,fs);
Plot_State_Colors(labels_frames,signif_frames,ColorState,OriginalExperiment,THR,fs,CoAc,Indexes);
Label_Condition_Raster(EnsembleName,R_Condition,fs);   % Labels
%   Sorted ***************************************************************
Plot_Raster_Ensembles(OriginalExperiment,Indexes(New_Order_Clustering),5,fs);
Plot_State_Colors(labels_frames,signif_frames,ColorState,OriginalExperiment,THR,fs,CoAc,Indexes(New_Order_Clustering));
Label_Condition_Raster(EnsembleName,R_Condition,fs);   % Labels
% Ensemble Transitions
% [ensemble_index_total]=Ensembles_Transitions(fs,labels_frames,signif_frames,ColorState,1);
% XY_cluster=XY_selectedClean(New_Order_Clustering,:);    