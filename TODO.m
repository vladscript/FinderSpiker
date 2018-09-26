%% NOTES##################################################################
% Manual Mode->Necessary to know the statistics power
% of the automatic method by dividing in -+ and -- (false+ & false-)
% Old Version Manual Mode:
%   Manual_Driver_Raster_Magic.m
%% FIXED  RECENTLY @ GIT
% BUG in Merge Selector: color doesn't update from contrast editor [SOLVED]
% Set color in DYE image AS WELL [DONE]
% NAVIGATE Making Zoom per EACH Coordinate [DONE]
% SAVE RoA Features at Select_Raster_for_NN

%% Bugs & New Functions NOW

% Show boxplots of Raster Features

% In the GUI of Ca++ Transients: update sparse signal->CLEANER when modify
%           lambda

% FEATURES of the MERGED NERUONS !!!
% Rate of ACtivity (frequentist probability)

% Add button to save Zoom image (MERGED MAGIC)
% Save Points SELECTED-> add to file .mat


% Add Highlight Neuron Using Mouse at Plot_Raster
% and other colors in the MERGE script : MAGENTA
% Threshold to get NETWORK !!!!!!!!
% Automatize Clustering
% Inspection of Each ROI...
% Detect when its empty detected or undetected at :
%   Undetected_Visual_Inspection


%% FUTURE **********************************
% Figure: reason whi mean(ROI) withput distortion
% Load Raw FLuorescenc vs F_0 distortion
% Analyze Rejects Ones Anyway to infer Artifacts
% Kalman Filtering at SNR and lambdas pdf's: for optimal threshold
% Automatize MERGE SELECTOR
% CLUSTERING STUFF ***********************
% Threshold: prior numbercoactivyt:
% [THCOAC]=mctest(R,'modes')
% Clustering
% [Ensembles,N_Ensembles,MethodsClustering,ThEffective,
% label_frames,
% signifi_frames]=ensemble_clusterin(R,THCOAC);

% Setup Intel/Info .mat File-> Default User Direcotry to save info
% Setup Script
% Detrending Issues
% High pass & Low Pass Filters :
%     > before detrending
%     > after detrending
%% STEPS ******************************************
%PROCESSING
% RUN >>Raster_Magic_Better
% RUN >>Detected_Visual_Inspection
% RUN >>Undetected_Visual_Inspection
% RUN >>PLot_and_Save


% RASTER SELECTION
% ACTUAL MODE: @ Original Coordiantes Order
% >>[RASTER_Selected_Clean,XY_selected,R_Condition,Onsets]= Select_Raster_for_NN(fs,Raster_Condition,XY,Names_Conditions,Experiment);

% >>Merge_Finder_Magic

% >>R=RASTER_Selected_Clean';
% >>R_CONDTION1=R_Condition{1}';
% ...
% >>R_CONDTIONi=R_Condition{i}';

% LOAD RASTER FEATURES
% >> Raster_Features_Display

% CLUSTERING
% 'Got to Dir' :
% C:\Users\Vladimir\Documents\Doctorado\Software\GetTransitum\Calcium Imaging Signal Processing\NeuralNetworks
% >>NeuralNetwork %-> GUI mode-> Clustering Analysis
% Back to 'FinderSPiker Dir '-
% >> Ensemble_Sorting

% >>Plot_Ensembles_Experiment(R_Condition,EnsembleName,Ensembled_Labels,Ensemble_Threshold,UniRMutiE,ColorState,fs,[]);
% >>Plot_Hebbian_Paths(R_Condition,Ensemble_Threshold,Ensembled_Labels,Names_Conditions,ColorState,fs);
% >>[Features_Ensemble,Features_Condition]=get_ensembles_features(R_Condition,Ensemble_Threshold,Ensembled_Labels,fs);
% >>save_features_ensembles(Experiment,Names_Conditions,Features_Ensemble,Features_Condition)


% COLOCALIZATION OF MARKED CELLS
% >>Merge_Finder_Magic




%% Legacy Code:
%   Manual_Driver_Raster_Magic_Ultimate (save as private mine)