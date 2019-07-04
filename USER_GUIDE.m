%% STEPS GUIDE ############################################################
% 1. SIGNAL PROCESSING: Detect Calcium Transients Events**********************
% >>Finder_Spiker_Calcium
% >>Detected_Visual_Inspection
% >>Undetected_Visual_Inspection
% >>Save_and_Plot
% REVIEW DATA:
% >>Plot_Experiment

% 2. RASTER SELECTION*********************************************************
% ACTUAL MODE: @ Original Coordiantes Order
% >>Select_Rasters
% % Save CSV Raster-Features of Merged & NO-Merged Cells:
% >>Select_Raster_PositiveCells;
% >>Select_Raster_NegativeCells
% 
% 2.1 CHECK RASTER's Selection DURATIONs: *************************************
% >>RasterDurations=get_raster_durations(Onsets,R_Condition,fs);
% 
% 2.2 TOTAL NETWORK (without Ensembles) ***************************************
% Save Links Features Without neither Ensembling nor Thresholding
% >>Get_Total_Network
% Show Boxplots and Make Table of Features: ONLY TOTAL or POSITIVE(!)
% >>Links_Features_Display
% 
% 2.3 RETRIEVE RASTER for ANALYSIS ********************************************
% >>R=RASTER_Selected_Clean'; % ALL CONDITIONS CONCATeNATED
% >>R_CONDITION1=R_Condition{1}; % Cells x Frames (dim)
% ...
% >>R_CONDITIONi=R_Condition{i};

% 3. CLUSTERING NEURONAL ENSEMBLES *******************************************
% AUTOMATICAL MAGIC
% >>R_CONDITIONi_Analysis=get_bayes_ensembles(R_CONDTIONi);
% HANDCRAFTED GUI (by JP): 
% >>NeuralNetwork 

% 3.1 DISPLAY AND SAVE RESULTS OF ENSEMBLES DISPLAY AND SAVE (GUI) ************
% Neural ensemble and Functional Network Features Extraction
% >> Ensemble_Sorting

% 3.2 PLOT ENSEMBLES FAST *****************************************************
% >> ImageEnsembles(R_ConditionNamej_Analysis);     % without  Hebbian Sequences
% >> ImageEnsembles(R_ConditionNamej_Analysis,1);   % with Hebbian Sequences

% 4. COLOCALIZATION OF MARKED CELLS ******************************************
% % Previously LOAD MAT FILE (?)
% >>Merge_Finder_Magic
%       It gets outputs: R_merged,R_nomerged,MetaDataColocaliation
% % Check Raster plots:
% >>Plot_Merged_NotMerged

% 5. RETRIEVE ORIGINAL SIGNALS from RASTER SELECTION ************************
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

% 6. ACCUMULATE FEATURES FROM SEVERAL EXPERIMENTS ****************************
% Choose One-by-One .mat Files-> Save .mat Files:
% >>Accumulate_Raster_Distances; % Save at CSV Files
% >>Accumulate_RoA_IEI_ED
% >>Accumulate_Ensembles_RoEn_IEnI_EnD
% >>Accumulate_Simm_Matrix

% 7. RASTER/ENSEMBLE/NETWORK FEATURES AND TABLES MAKER ***********************
% 
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
% >>Merge_Feature_Datasets          % Make DATASET of ONE Kind MULTI EXPS
%  Get ALL Features Dataset: Raster,Ensembles,Networks
% >>Merge_Datasets                  % Make DATSET of the 3 Kinds


% 8. DATA FEATURE EXPLORING: After Merging Datasets or Single Type Datasets.
% >>Feature_Explorer         *UNDER CONSTRUCTION*

% 9. MACHINE LEARNING: choose a Dataset:
% >>Features_Datasets_NBC; % *UNDER CONSTRUCTION*


%% END ####################################################################