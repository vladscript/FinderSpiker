%% PERMANTENT NOTES # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
% Manual Mode is divideed in two parts->
% Necessary to know the statistical power of the Automatic Mode:
% of the automatic method by dividing in -+ and -- (false+ & false-)

%% FIXED  READY TO GO @ GIT

% >>Accumulate_Raster_Distances make ti to recognize Merge or NotMerged
% Directories to save Network Stuff

%% Bugs & New Functions NOW

% Statistics: automate: OK

% Features From Get_Total_Network?
% Links Features Without->Merge THEM

% Adjust Contrast and Sliding Cells in Merge_Finder

% Dyskinesia Data: Get_Total_Network

% RETRIEVE SIGNALS & RECONSTRUCT VIDEO
% retreive of Original Signals, coordinates, etc:
% Re make clean video

% delete >Plot_Raster_V.m; gen_feat_table_merged.m; get_merged_coordinates;
% NN2Gephi.m; Raster2fNet.m

% Make PCA of RASTER: 
% denoise raster: get most variance PCs and rebuild raster
% Add Colocalizer Filters: TdTomato, Yellow, Others
%  Update to MATLAB 2019a !!!!! UNAM License
% Test Mike's Clustering Algorithm 4 CoActivity (transposed matrix)and save as weel

%%% MAKE ALGORITHMIA
% Add button to save Zoom image (MERGED MAGIC)
% Save Selected Points SELECTED-> add to file .mat
% Add Highlight Neuron Using Mouse at Plot_Raster
% and other colors in the MERGE script : MAGENTA
% Inspection for Each ROI...
% Make RGB sum at Functional Network
% driver ghost Issuee(?)
% MAKE NETWORK highlight special neuron population

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
% Save Links Features Without neither Ensembling nor Thresholding
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
% >> ImageEnsembles(R_ConditionNamej_Analysis);     % without  Hebbian Sequences
% >> ImageEnsembles(R_ConditionNamej_Analysis,1);   % with Hebbian Sequences

% COLOCALIZATION OF MARKED CELLS ******************************************
% % % Previously LOAD MAT FILE (?)
% 
% >>Merge_Finder_Magic
%       It gets outputs: R_merged,R_nomerged,MetaDataColocaliation
% 
% % % Save CSV Raster-Features of Merged & NO-Merged Cells:
% 
% >>Select_Raster_for_NN(fs,R_merged,XY,Names_Conditions,Experiment);
% >>Select_Raster_for_NN(fs,R_nomerged,XY,Names_Conditions,Experiment);
% 
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

% ACCUMULATE FEATURES FROM SEVERAL EXPERIMENTS ****************************
% Choose One-by-One .mat Files-> Save .mat Files:
% >>Accumulate_Raster_Distances; % Save at CSV Files
% >>Accumulate_RoA_IEI_ED
% >>Accumulate_Ensembles_RoEn_IEnI_EnD
% >>Accumulate_Simm_Matrix

% RASTER/ENSEMBLE/NETWORK FEATURES AND TABLES MAKER ***********************
% RASTER FEATURES
% >>Raster_Features_Display                 % make CSV table
% ENSEMBLES GENERAL FEATURES
% >>Ensembles_Features_Display              % make CSV table
% ENSEMBLES DETAILED FEATURES 
% >>Ensembles_Features_Detailed_Display     % make CSV table
% MAKE DATASETS FROM GEPHI NETWORK FEATURES
% >>Get_Gephi_Data;                 % Make .mat Dataset from Gephi Files
% >>Make_Statistics_Gephi_Features; % Make csv table of SINGLE Feature
% >>Concatenate_NetFeats            % Make table of MULTIPLE Features
% 
% Script to Merged them and make a DATASET for Machine Learning:
% >>Merge_Feature_Datasets          % Make DATASET of ONE Kind
%  Get ALL Features Dataset: Raster,Ensembles,Networks
% >>Merge_Datasets                  % Make DATSET of the 3 Kinds
% 

% DATA FEATURE EXPLORING: After Merging Datasets or Single Type Datasets.
% >>Feature_Explorer         *UNDER CONSTRUCTION*

% MACHINE LEARNING: choose a Dataset:
% >>Features_Datasets_NBC; % *UNDER CONSTRUCTION*


%% END ####################################################################