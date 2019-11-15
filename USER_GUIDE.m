%% SOURCE
% $FinderSpiker$
% 
% <https://github.com/vladscript/FinderSpiker.git>
% 
%%#########################################################################
% ######################### *STEPS GUIDE* #################################
% #########################################################################
% 

%% 0. SETTTINGS
% 
% * Deconvolution Parameters:   >>edit Load_Default_Values_SP
% * Clustering Parameters:      >>edit Load_Default_Clustering
% * Colors:                     >>edit SetColorMap.m
% 
% 
%% 1. SIGNAL PROCESSING: Detect Calcium Transients Events
% Getting Activation Matrix (raster) from Calium Imaging Data for each Cell
% DATA:                 VIDEOS, Sampling Frequency & Set of Coordinates
% 
% * >>Finder_Spiker_Calcium
% * >>Detected_Visual_Inspection
% * >>Undetected_Visual_Inspection
% * >>Save_and_Plot
% 
% REVIEW DATA:
% 
% * >>Plot_Experiment
% 
%% 2. RASTER SELECTION
% ACTUAL MODE: @ Original Coordiantes Order
% 
% * >>Select_Rasters
% 
% Save CSV Raster-Features of Merged & NO-Merged Cells (Step 4):
% 
% * >>Select_Raster_PositiveCells;
% * >>Select_Raster_NegativeCells
% 
% 2.1 CHECK RASTER's Selection DURATIONs: *********************************
% 
% * >>RasterDurations=get_raster_durations(Onsets,R_Condition,fs);
% 
% 2.2 TOTAL NETWORK (without Step 3) **************************************
% 
% Save Links Features Without neither Ensembling nor Thresholding
% 
% * >>Get_Total_Network
% 
% Show Boxplots and Make Table of Features: ONLY TOTAL or POSITIVE(!)
% 
% * >>Links_Features_Display
% 
% 2.3 RETRIEVE RASTER for ANALYSIS ****************************************
% 
% * >>R=RASTER_Selected_Clean'; % ALL CONDITIONS CONCATeNATED
% * >>R_CONDITION1=R_Condition{1}; % Cells x Frames (dim)
% ...
% * >>R_CONDITIONi=R_Condition{i};
% 2.4 Plot Sorted by Rate of Activiyt
% 
% * >>[IndexSorteActivity]=plot_activityrate(R_Condition,Names_Conditions,fs);
% 
%% 3. CLUSTERING NEURONAL ENSEMBLES
% AUTOMATICAL & MANUAL MAGIC
% 
% * >>R_CONDITIONi_Analysis=get_bayes_ensembles(R_CONDTIONi);
% 
% * >>R_CONDITIONi_Analysis=get_ensembles(R_CONDTIONi,CAGthreshold,Nensambles);
% 
% HANDCRAFTED GUI (by JP): 
% 
% * >>NeuralNetwork 
% 
% 3.1 DISPLAY AND SAVE RESULTS OF ENSEMBLES DISPLAY AND SAVE (GUI) ********
% Neural ensemble and Functional Network Features Extraction
% 
% * >> Ensemble_Sorting
% 
% 3.2 PLOT ENSEMBLES FAST *************************************************
% 
% * >> ImageEnsembles(R_ConditionNamej_Analysis,1);   %with Hebbian Sequences
% * >> ImageEnsembles(R_ConditionNamej_Analysis);     %without  
% 
%% 4. COLOCALIZATION OF MARKED CELLS 
% % Previously LOAD MAT FILE (?)
% 
% * >>Merge_Finder_Magic
% 
%       It gets outputs: R_merged,R_nomerged,MetaDataColocaliation
% % Check Raster plots:
% 
% * >>Plot_Merged_NotMerged
% 
%% 5. RETRIEVE ORIGINAL SIGNALS from RASTER SELECTION 
% 5.1 Display Signals:
% 
% 0) Get Merged Coordinates (IF SO)
% 
% * >> XY_merged=XY_selected(MetaDataColocaliation.PositiveCells,:);
% 
% 1) Plot Raster (*without sorting*) from:
% 
% * >>[Rsel,IndexSorted]=Retrieve_Selected_Signal(Onsets,R_Condition,RASTER,XY_subset,XY);
% * >>Rjunto=[Rsel{1},Rsel{...},Rsel{NConditions}]; 
% * >>Plot_Raster_Ensembles(Rjunto,fs)
% * >>Label_Condition_Raster(Names_Conditions,R_Condition,fs); 
% 
% 2) Find Cell Signal of Interest: Ci
% 
% * >>[XS,IndexSorted]=Retrieve_Selected_Signal(Onsets,R_Condition,SIGNALSclean,XY_subset,XY);
% * >>figure; plot(XS{c}(Ci,:))
% 
% 5.2 Clean Video Output:
% 
% * >>Video_Cleaner;    % Edit @ Script
% 
%% N-EXPERIMENTS RESULTS 
% 
%% 6. ACCUMULATE FEATURES FROM SEVERAL EXPERIMENTS
% 
%   Choose One-by-One .mat Files-> Save .mat Files:
% 
% * >>Accumulate_Raster_Distances; % Save at CSV Files
% * >>Accumulate_RoA_IEI_ED
% * >>Accumulate_Ensembles_RoEn_IEnI_EnD
% * >>Accumulate_Simm_Matrix
% 
%% 7. DATASETS MAKER 
% 
% Make Tables for SINGLE EXPERIMENTS
% 
% 7.1 RASTER FEATURES
% 
% * >>Raster_Features_Display                 % make CSV table
% 
%                                           % Make Special Directories
% 7.2 ENSEMBLES GENERAL FEATURES
% 
% * >>Ensembles_Features_Display              % make CSV table
% 
%                                           % Make Special Directories
% 7.3 ENSEMBLES DETAILED FEATURES 
% 
% * >>Ensembles_Features_Detailed_Display     % make CSV table
% 
%                                           % Make Special Directories
% 7.4 LINKS FEATURES
% 
% * >>Links_Features_Display                  % ?
% 
%% 8 MAKE DATASETS FROM GEPHI NETWORK FEATURES
% 
% 7.9  Export 'Workspaces' from Gephi     Data could be in  '\NetWorks-CSV'
% 
%                                         Using Gephi 0.9.1
% 8.1  Make Network Features Dataset
% 
% * >>Get_Gephi_Data;      % Make .mat Dataset from CSV Gephi Workspaces
% 
%                        Save .mat File @ FinderSpiker/DatabaseFolder
% 
% * >>Make_Statistics_Gephi_Features; 
% 
%                        Read .mat File @ FinderSpiker/DatabaseFolder
%                        Show RAINPLOT for each Feature     
%                        Create Network_FEATURE_Dataset
% 
% * >>Concatenate_NetFeats            
% 
%                       Read Network_FEATURE_Dataset
%                       Concatenate table of MULTIPLE Features in one
%                       Create NetFEATS_Dataset @FinderSpiker/DatabaseFolder
% 
%% 9  STACK DIVERS DATASETS
% 
% % Make DATASET for Stack Several tables from Divers Experiments:
% 
% * >>Merge_Feature_Datasets          
% 
%   -Raster Activity:       Table
%   -General Ensembles:     Tables 
%   -Detailed Ensembles:    Tables 
%   -Network Features:      NetFEATS 
% 
% E.G. DIVERS set of Experiments,i.e, CTRL, PARKINSON, DYSKINESIA, ...
% 
%% 10 DATA FEATURE EXPLORER: Merging Feature Datasets 
% 
% 10.1  Get ALL Features Dataset: Raster,Ensembles,Networks
% 
% * >>Merge_Datasets           % Make DATASET of the 3 Kinds in 1 file
%                               (not needed for Step 10.2)
% 
% 10.2  Dysplay All Features RAINPLOTS & Statistics
% 
% * >>Feature_Explorer         *UNDER CONSTRUCTION*
% 
% 
%% 11 MACHINE LEARNING: choose a Dataset:
% 
% * >>Features_Datasets_NBC; % *UNDER CONSTRUCTION*
% 

%% END ####################################################################