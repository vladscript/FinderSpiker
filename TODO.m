%% PERMANENT NOTES # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
% Manual Mode is divided in two parts->Automatic and Standard:
% [ Necessary to know the statistical power of the Automatic Mode,
% Getting insohgts about -+ and -- (false+ & false-)]
% Transients are calculated from denoised signal (wavelet analysis)
% Rate of Activity is from Actual and Previou(s) conditions only
% Get Only Ensembles that last at least 1/2 second
% Hebbian sequence holds if CAG level stands still
%           TypeCycles(1)->Simple
%           TypeCycles(2)->Closed
%           TypeCycles(3)->Open
% Network Link Weight Threshold: Minimum RoA of Neuron to Belong any ensemble
%% FIXED  READY TO GO @ GIT
%   Ensemble Intervals Fix
% READ ROIs from ImageJ -----------------

%% Bugs & New Functions NOW


% MAKE Plot_Simm_Matrix.m from labeled_frames vector

% Make Experimental Datasets Folder for 'Merge_Feature_Datasets'

% Zero-Cluster issue at 'get_ensembles'

% Dissapeared Ensembles by labeling Sequence: see 20190930-1 
% R_DAdepleted_Analysis=get_ensembles(R_DAdepleted',10,3);

% Statistics: automate
% CAG threshold (Mike's Script - MonteCarlo)
% Recover Coordinates of Ensembles and Display as scatter plots

% INDEX of USEFUL FUNCTIONS

% Update to MATLAB 2019a ! UNAM License[other versions can be downloaded]

% Toy example for clustering: time correlaition vs coactivity:
%  to seen syn fires is irrelevante
% coactivity: measurebale (Published) features
% correlation: to publish....

% Make PCA of RASTER: 
% denoise raster: get most variance PCs and rebuild raster


% BEFORE TO RELEASE *******************************************************
% DELETE >Plot_Raster_V.m; gen_feat_table_merged.m; get_merged_coordinates;
% NN2Gephi.m; Raster2fNet.m; get_iti_pdf(?); interval_duration_events(?);
% Plot_Calcium_Transients.m; get_ensembles_manual.m
% Plot_Raster_Results; get_ensemble_intervals; Plot_Ensembles;

%% FUTURE *****************************************************************
%%% MAKE ALGORITHMIA
% Export Rasters @ CSV format
% Adjust Contrast and Sliding Cells in Merge_Finder
% Normalize Amplitude @ Detrended Signal
% Add button to save Zoom image (MERGED MAGIC)
% Save Selected Points SELECTED-> add to file .mat
% Add Highlight Neuron Using Mouse at Plot_Raster
% and other colors in the MERGE script : MAGENTA
% Inspection for Each ROI...
% driver ghost Issuee(?)->Response Size
% MAKE NETWORK highlight special neuron population
% 
% Look at Line Equation: 0:N-1 or 1:N   [*]
% Figure: reason whi mean(ROI) without distortion
% Load Raw FLuorescenc vs F_0 distortion
% Analyze Rejects Ones Anyway to infer Artifacts
% Processing Times/Detections/etc from log files
% Automatize MERGE SELECTOR
% Setup Intel/Info .mat File-> Default User Directory to save info
% Setup Script: deconvolution parameters
% Check at Signals with Huge Valley (synaptic like)
% 
% Merge tables-> if conditions ain't the same
% Add Colocalizer Filters: TdTomato, Yellow, Others
%% ########################################################################