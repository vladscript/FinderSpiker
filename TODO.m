%% PERMANENT NOTES # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
% Manual Mode is divided in two parts->Automatic and Standard:
% [ Necessary to know the statistical power of the Automatic Mode,
% Getting insohgts about -+ and -- (false+ & false-)]
% Transients are calculated from denoised signal (wavelet analysis)
% Rate of Activity Network Features is from Actual and Previou(s) conditions only
% Get Only Ensembles that last at least 1/2 second
% Hebbian sequence holds if CAG level stands still
%           TypeCycles(1)->Simple
%           TypeCycles(2)->Closed
%           TypeCycles(3)->Open
% Network Link Weight Threshold: Minimum RoA of Neuron to Belong any ensemble
% Use Nodes with Degree>0 in C_i or peviously Degree>0 in C_i-1
% 
%% URGENT
% 
% EDIT to load data from functions: Display_NetworkPDFs !!!!!!!!!!!!
% Datos de neurodesarrollo UPDATE
% 
%% READY TO GO:

% Fucntionality MAIN updates:
% HELP GUI for Calcium Magic
% Segmentation of CSV Files
% Display_NetworkPDFs Histograms with Lines instead of bars
% Hebbian Times in Feature_Condition strcut
% pairedraincloud UPDATED to read Cell arrays
% Network PDFs updates @ Display_NetworkPDFs
% Classification_Regions: KernelScale:fixed and OVA balarcolor normalized
% Included IndexSorting for Plot by activity rate
)
% 
%% TO DO:
%  INITIAL MENU: ImPatch Data or CSV Data
% READ AND ANALYSE CSV Data
% SVM paramters @ svmmulticlass
% Import Graphs in R (igraph->
% Intro interface @ Classification_Regions [PCA->SVM, MC & OVA stuff]
% Load modesl and evaluate new Data for ML:
%  It requires make new table of new data

% Check Merge_Feature_Datasets-> maybe obsolete

% Publish USER_GUIDE
% Feature_Explorer with Statistics and cloud plots (?)

% Plot_Experiment: Load .mat File if Wrokspace "Empty"
% Export Rasters to CSVs
% 
% USE: https://umap-learn.readthedocs.io/en/latest/basic_usage.html
%% Bugs & New Functions NOW

% MAKE Plot_Simm_Matrix.m from labeled_frames vector

% Make Experimental Datasets Folder for 'Merge_Feature_Datasets'

% Zero-Cluster issue at 'get_ensembles' (?)

% Dissapeared Ensembles by labeling Sequence: see 20190930-1 
% R_DAdepleted_Analysis=get_ensembles(R_DAdepleted',10,3);

% Recover Coordinates of Ensembles and Display as scatter plots

% INDEX of USEFUL FUNCTIONS

% Toy example for clustering: time correlaition vs coactivity:
%  to seen syn fires is irrelevante

% Make PCA of RASTER: 
% denoise raster: get most variance PCs and rebuild raster


% BEFORE TO RELEASE *******************************************************
% DEPRECATE> ;
% get_iti_pdf(?); interval_duration_events(?);
% 

%% FUTURE *****************************************************************
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