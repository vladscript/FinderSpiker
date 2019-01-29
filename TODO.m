%% PERMANTEN NOTES # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
% Manual Mode is divideed in two parts->
% Necessary to know the statistical power of the Automatic Mode:
% of the automatic method by dividing in -+ and -- (false+ & false-)
% Old Version Manual Mode:
%   Manual_Driver_Raster_Magic.m (becoming unnecessary)
%% FIXED  READY TO GO @ GIT

%% Bugs & New Functions NOW

% Look for Bugs:
% Ignore Single-Frame Ensembles-> OK (zero vector)
% Bugs at Rasters with Very LOW ACTIVITY Exp?
% Save .mat File even when if it was analyzed at NeuralNetwork GUI

% Compare SIMMILAR and CONFUSING (confmat) ENSEMBLES to MERGE them
% Merge the Most inter-Cinfusing if its distance it's close enough

% SIGNAL PROCESSING RELATED
%   Somee offset at detrending algorithm 
%   Spurious Drivers
%       Lone Drivers: check clean signal's samples around if they're above noise
%   Check at Signals with Huge Valley (synaptic like)

% retreive of Original Signals, coordinates, etc:

%%% MAKE ALGORITHMIA

% Modify Zoom at Windows in the script: Select_Raster_for_NN.m

% delete Plot_Raster_V.m

% Add button to save Zoom image (MERGED MAGIC)
% Save Selected Points SELECTED-> add to file .mat
% Add Highlight Neuron Using Mouse at Plot_Raster
% and other colors in the MERGE script : MAGENTA

% Threshold to get NETWORK !!!!!!!!

% Inspection of Each ROI...
% Detect when its empty detected or undetected at :
%   Undetected_Visual_Inspection

% CLUSTERING STUFF ***********************
% Threshold: Given by Best Classifiaction Validation


%% FUTURE **********************************
% Figure: reason whi mean(ROI) withput distortion
% Load Raw FLuorescenc vs F_0 distortion
% Analyze Rejects Ones Anyway to infer Artifacts
% Automatize MERGE SELECTOR
% Setup Intel/Info .mat File-> Default User Directory to save info
% Setup Script: deconvolution parameters

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
% MANUAL
% >>R_COND_Analysis = raster_cluster(R_COND,CAG_TH,Nensmbles,'hamming');

% DISPLAY AND SAVE RESULTS OF ENSEMBLES DISPLAY AND SAVE (GUI)
% Neural ensemble and Functional Network Features Extraction
% >> Ensemble_Sorting

% PLOT FAST
% >> ImageEnsembles(R_ConditionNamej_Analysis);

% COLOCALIZATION OF MARKED CELLS
% >>Merge_Finder_Magic

% INPSECTION to RETRIEVE ORIGINAL SIGNALS from RASTER SELECTION
% 1) Plot Raster (without sorting) from:
% >>[Rsel,IndexSorted]=Retrieve_Selected_Signal(Onsets,R_Condition,RASTER,XY_merged,XY);
% Each Rsel{n} is the selected raster
% >>Rjunto=[Rsel{1},Rsel{2},Rsel{3}]; Plot_Raster_Ensembles(Rjunto);
% Find Cell of Interest: Ci
% >>[XS,IndexSorted]=Retrieve_Selected_Signal(Onsets,R_Condition,SIGNALSclean,XY_merged,XY);
% >>figure; plot(XS{c}(Ci,:))

% N-EXPERIMENTS RESULTS ******

% LOAD RASTER FEATURES FOR A SET OF EXPERIMENTS
% >> Raster_Features_Display

% ENSEMBLES FEATURES BOXPLOTS (*)
% >>Ensembles_Features_Display

% Cummulative Features for differente Features:
% such as RoA of several cells and experiments

%% Legacy Code:
%   Manual_Driver_Raster_Magic_Ultimate (saved as private)