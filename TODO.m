%% PERMANTENT NOTES # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
% Manual Mode is divideed in two parts->
% Necessary to know the statistical power of the Automatic Mode:
% of the automatic method by dividing in -+ and -- (false+ & false-)
% Old Version Manual Mode:
%   Manual_Driver_Raster_Magic.m (becoming unnecessary)
%% FIXED  READY TO GO @ GIT
% Stats from Any Netwrok Feature: done!

%% Bugs & New Functions NOW

% STATS visualization!
% Statistics: automate
% Hamming Distance among active Neurons of Different Conditions

% ANALYZE DYSKINESIA DATA

% MAKE NETWORK highlight special neuron population


% Test Mike's Clsutering Algorithm 4 CoActivity (transposed matrix)and save as weel

% RETRIEVE SIGNALS & RECONSTRUCT VIDEO
% retreive of Original Signals, coordinates, etc:
% Re make clean video

%%% MAKE ALGORITHMIA

% delete >Plot_Raster_V.m; gen_feat_table_merged.m; get_merged_coordinates;

% Make PCA of RASTER: 
% denoise raster: get most variance PCs and rebuild raster

% Add button to save Zoom image (MERGED MAGIC)
% Save Selected Points SELECTED-> add to file .mat
% Add Highlight Neuron Using Mouse at Plot_Raster
% and other colors in the MERGE script : MAGENTA
% Inspection for Each ROI...
% Make RGB sum at Functional Network
% driver ghost Issuee(?)

%% FUTURE *****************************************************************
% Driver Issues
% Look at Line Equation: 0:N-1 or 1:N   [*]
% Figure: reason whi mean(ROI) without distortion
% Load Raw FLuorescenc vs F_0 distortion
% Analyze Rejects Ones Anyway to infer Artifacts
% Processing Times/Detections/etc from log files
% Automatize MERGE SELECTOR
% Setup Intel/Info .mat File-> Default User Directory to save info
% Setup Script: deconvolution parameters
% Check at Signals with Huge Valley (synaptic like)

%% STEPS GUIDE ############################################################

% SIGNAL PROCESSING: Detect Calcium Transients Events**********************
% >>Finder_Spiker_Calcium
% >>Detected_Visual_Inspection
% >>Undetected_Visual_Inspection
% >>Save_and_Plot
% REVIEW DATA:
% >>Plot_Experiment

% RASTER SELECTION*********************************************************
% ACTUAL MODE: @ Original Coordiantes Order
% >>Select_Rasters

% CHECK RASTER's Selection DURATIONs: *************************************
% >>RasterDurations=get_raster_durations(Onsets,R_Condition,fs);

% TOTAL NETWORK (without Ensembles) ***************************************
% Save Links Features without Thresholding
% >>Get_Total_Network

% RETRIEVE RASTER for ANALYSIS ********************************************
% >>R=RASTER_Selected_Clean'; % ALL CONDITIONS CONCATeNATED
% >>R_CONDITION1=R_Condition{1}; % Cells x Frames (dim)
% ...
% >>R_CONDITIONi=R_Condition{i};

% CLUSTERING NEURONAL ENSEMBLES *******************************************
% AUTOMATIC
% >>R_CONDITIONi_Analysis=get_bayes_ensembles(R_CONDTIONi);
% MANUAL GUI (by JP): 
% >>NeuralNetwork 

% DISPLAY AND SAVE RESULTS OF ENSEMBLES DISPLAY AND SAVE (GUI) ************
% Neural ensemble and Functional Network Features Extraction
% >> Ensemble_Sorting

% PLOT ENSEMBLES FAST *****************************************************
% >> ImageEnsembles(R_ConditionNamej_Analysis); % without  Hebbian Sequences
% >> ImageEnsembles(R_ConditionNamej_Analysis,1); % with Hebbian Sequences

% COLOCALIZATION OF MARKED CELLS ******************************************
% % % Previously LOAD MAT FILE (?)
% >>Merge_Finder_Magic
%       It gets outputs: R_merged,R_nomerged,MetaDataColocaliation
% Save CSV Raster-Features of Merged & NO-Merged Cells:
% >>Select_Raster_for_NN(fs,R_merged,XY,Names_Conditions,Experiment);
% >>Select_Raster_for_NN(fs,R_nomerged,XY,Names_Conditions,Experiment);
% Check Raster plots:
% >>Plot_Merged_NotMerged

% RETRIEVE ORIGINAL SIGNALS from RASTER SELECTION ************************
% 0) Get Merged Coordinates (IF SO)
% >> XY_merged=XY_selected(MetaDataColocaliation.PositiveCells,:);
% 1) Plot Raster (*without sorting*) from:
% >>[Rsel,IndexSorted]=Retrieve_Selected_Signal(Onsets,R_Condition,RASTER,XY_subset,XY);
% >>Rjunto=[Rsel{1},Rsel{...},Rsel{NConditions}]; 
% >>Plot_Raster_Ensembles(Rjunto,fs)
% >>Label_Condition_Raster(Names_Conditions,R_Condition,fs); 
% 2) Find Cell Signal of Interest: Ci
% >>[XS,IndexSorted]=Retrieve_Selected_Signal(Onsets,R_Condition,SIGNALSclean,XY_subset,XY);
% >>figure; plot(XS{c}(Ci,:))

% N-EXPERIMENTS RESULTS ###################################################

% LOAD & GET FEATURES FOR A SET OF EXPERIMENTS*****************************
% These scripts save feature tables 
% for Machine Learning or Statistical Analysis

% RASTER FEATURES AND TABLES MAKER
% Choose All CSV files at once
% >>Raster_Features_Display
% ENSEMBLES GENERAL FEATURES 
% >>Ensembles_Features_Display
% ENSEMBLES DETAILED FEATURES 
% >>Ensembles_Features_Detailed_Display
% Script to Merged them and make a DATASET for Machine Learning:
% >>Merge_Feature_Datasets

% MAKE DATASETS FROM GEPHI NETWORK FEATURES
% >>Get_Gephi_Data;
% >>Make_Statistics_Gephi_Features;

% MACHINE LEARNING: choose a Dataset:
% >>Features_Datasets_NBC;

% DATA EXPLORING
% >>Feature_Explorer

% ACCUMULATE FEATURES FROM SEVERAL EXPERIMENTS ****************************
% Choose One-by-One .mat Files-> Save .mat Files:
% >>Accumulate_RoA_IEI_ED
% >>Accumulate_Ensembles_RoEn_IEnI_EnD
% >>Accumulate_Simm_Matrix

% Plot_Accumulate_CDF; % To plot Results
%% END ####################################################################