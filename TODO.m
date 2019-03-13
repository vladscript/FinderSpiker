%% PERMANTENT NOTES # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
% Manual Mode is divideed in two parts->
% Necessary to know the statistical power of the Automatic Mode:
% of the automatic method by dividing in -+ and -- (false+ & false-)
% Old Version Manual Mode:
%   Manual_Driver_Raster_Magic.m (becoming unnecessary)

% Features  of Experiments:  
% >Rasters,             MATLAB script
% >Ensembles general    MATLAB script
% >Ensembles details    MATLAB script
% >Networks             from Gephi

%% FIXED  READY TO GO @ GIT
% Big Bug @ Detrending Algortithm

%% Bugs & New Functions NOW

% Get best subset of features that best classify vs PCA

% All-Features->PCA (dim red)-> SVM: not really good

% Statistics

% Test Mike's Clsutering Algorithm 4 CoActivity (transposed matrix)and save as weel


% RETRIEVE SIGNALS & RECONSTRUCT VIDEO
% retreive of Original Signals, coordinates, etc:
% Re make clean video

%%% MAKE ALGORITHMIA

% delete >Plot_Raster_V.m >Plot_Merged_NotMerged.m

% Make PCA of RASTER: 
% denoise raster: get most variance PCs and rebuild raster

% Test Visualizer of CDF for (+)and (-) colocated cells

% Add button to save Zoom image (MERGED MAGIC)
% Save Selected Points SELECTED-> add to file .mat
% Add Highlight Neuron Using Mouse at Plot_Raster
% and other colors in the MERGE script : MAGENTA
% Inspection for Each ROI...


%% FUTURE **********************************
% Look at Line Equation: 0:N-1 or 1:N   [*]
% Figure: reason whi mean(ROI) withput distortion
% Load Raw FLuorescenc vs F_0 distortion
% Analyze Rejects Ones Anyway to infer Artifacts
% Processing Times/Detections/etc from log files
% Automatize MERGE SELECTOR
% Setup Intel/Info .mat File-> Default User Directory to save info
% Setup Script: deconvolution parameters
% Check at Signals with Huge Valley (synaptic like)

%% STEPS GUIDE *********************************************************
% SIGNAL PROCESSING: Detect Calcium Transients Events
% >>Finder_Spiker_Calcium
% >>Detected_Visual_Inspection
% >>Undetected_Visual_Inspection
% >>Save_and_Plot

% RASTER SELECTION
% ACTUAL MODE: @ Original Coordiantes Order
% >>Select_Rasters

% RETRIEVE RASTER for ANALYSIS
% >>R=RASTER_Selected_Clean'; % ALL CONDITIONS
% >>R_CONDITION1=R_Condition{1}; % Cells x Frames (dim)
% ...
% >>R_CONDITIONi=R_Condition{i};

% CLUSTERING NEURONAL ENSEMBLES
% AUTOMATIC
% >>R_CONDITIONi_Analysis=get_bayes_ensembles(R_CONDTIONi);
% MANUAL GUI (by JP): 
% >>NeuralNetwork 


% DISPLAY AND SAVE RESULTS OF ENSEMBLES DISPLAY AND SAVE (GUI)
% Neural ensemble and Functional Network Features Extraction
% >> Ensemble_Sorting

% PLOT ENSEMBLES FAST
% >> ImageEnsembles(R_ConditionNamej_Analysis); % without  Hebbian Sequences
% >> ImageEnsembles(R_Dyskinesia_Analysis,1); % with Hebbian Sequences

% COLOCALIZATION OF MARKED CELLS
% >>Merge_Finder_Magic
% Save CSV Features of Merged Cells:
% >>Select_Raster_for_NN(fs,R_merged,XY,Names_Conditions,Experiment);
% >>Select_Raster_for_NN(fs,R_nomerged,XY,Names_Conditions,Experiment);

% INSPECTION to RETRIEVE ORIGINAL SIGNALS from RASTER SELECTION
% 1) Plot Raster (without sorting) from:
% >>[Rsel,IndexSorted]=Retrieve_Selected_Signal(Onsets,R_Condition,RASTER,XY_merged,XY);
% Each Rsel{n} is the selected raster
% >>Rjunto=[Rsel{1},Rsel{2},Rsel{3}]; Plot_Raster_Ensembles(Rjunto);
% Find Cell of Interest: Ci
% >>[XS,IndexSorted]=Retrieve_Selected_Signal(Onsets,R_Condition,SIGNALSclean,XY_merged,XY);
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

% IGNORE EXPERIMENTS:
% >>Check_Data_Dyskinesia % (user defined)

% MACHINE LEARNING: choose a Dataset:
% >>Features_Datasets_NBC

% DATA EXPLORING
% >>Feature_Explorer

% ACCUMULATE FEATURES FROM SEVERAL EXPERIMENTS ****************************
% Choose One-by-One .mat Files
% >>Accumulate_RoA_IEI_ED
% >>Accumulate_Ensembles_RoEn_IEnI_EnD
% >>Accumulate_Simm_Matrix

%% END ####################################################################