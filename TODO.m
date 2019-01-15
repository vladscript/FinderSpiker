%% NOTES # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
% Manual Mode->Necessary to know the statistics power
% of the automatic method by dividing in -+ and -- (false+ & false-)
% Old Version Manual Mode:
%   Manual_Driver_Raster_Magic.m
%% FIXED  READY TO GO @ GIT
% Plot_Raster_V-> square size=4 (Prohibittedsize=3 or size=1)
%   make the dots_> squares

%% Bugs & New Functions NOW

% UNLABEL frames of CAG(n)=1 AND it belong to more than 1 ensemble
% PLOTTING COLORS: missing some fast alternated ensembles :




% Somee offset at detrending algorithm 

% Spurious Drivers
% Lone Drivers: check clean signal's samples around if they're above noise
% Check at Signals with Huge Valley (synaptic like)

% Modify Zoom at Windows in the script: Select_Raster_for_NN.m

% Bugs at Rasters with Very LOW ACTIVITY Exp?
% > Get Criteria for CAG_th & Nensembles->Partially Options A,B
% Get and Save CAG temporal Features: Autocorrelation, etc...
% IDENTIFY Features per CONDITION when Analyze Concatenated Raster
% (condition times)

% Save .mat File even when if it was analyzed at NeuralNetwork GUI

% Matrix of Dist: Hamming distance of neurons sets at each Condition
% Save Ensemble Features

% Add button to save Zoom image (MERGED MAGIC)
% Save Selected Points SELECTED-> add to file .mat
% Add Highlight Neuron Using Mouse at Plot_Raster
% and other colors in the MERGE script : MAGENTA
% Threshold to get NETWORK !!!!!!!!

% Inspection of Each ROI...
% Detect when its empty detected or undetected at :
%   Undetected_Visual_Inspection

% CLUSTERING STUFF ***********************
% Threshold: prior numbercoactivyt:
% [THCOAC]=mctest(R,'modes')


%% FUTURE **********************************
% Figure: reason whi mean(ROI) withput distortion
% Load Raw FLuorescenc vs F_0 distortion
% Analyze Rejects Ones Anyway to infer Artifacts
% Automatize MERGE SELECTOR
% Setup Intel/Info .mat File-> Default User Directory to save info
% Setup Script: deconvolution parameters

%% STEPS GUIDE *********************************************************
% PROCESSING
% >>Finder_Spiker_Calcium
% >>Detected_Visual_Inspection
% >>Undetected_Visual_Inspection
% >>Save_and_Plot

% RASTER SELECTION
% ACTUAL MODE: @ Original Coordiantes Order
% >>[RASTER_Selected_Clean,XY_selected,R_Condition,Onsets]= Select_Raster_for_NN(fs,Raster_Condition,XY,Names_Conditions,Experiment);

% RETRIEVE RASTER for ANALYSIS
% >>R=RASTER_Selected_Clean'; % ALL CONDITIONS
% >>R_CONDTION1=R_Condition{1}';
% ...
% >>R_CONDTIONi=R_Condition{i}';

% CLUSTERING ENSEMBLES
% For everi i conditoneach R_CONDITIONi:
% >>R_ConditionName1=R_Condition{1};
% AUTOMATIC
% >>R_ConditionName1_Analysis=get_bayes_ensembles(R_Dyskinesia);
% MANUAL
% Set CAG_TH,  Nensembles & SimMethod='hamming'
% >>Ri_Analysis = raster_cluster(R_CONDITIONi,CAG_TH,Nensembles,SimMethod); 

% DISPLAY AND SAVE RESULTS OF ENSEMBLES DISPLAY AND SAVE (GUI)
% >> Ensemble_Sorting

% PLOTSSS
% >>Plot_Ensembles_Experiment(R_Condition,EnsembleName,Ensembled_Labels,Ensemble_Threshold,UniRMutiE,ColorState,fs,[]);
% >>Plot_Hebbian_Paths(R_Condition,Ensemble_Threshold,Ensembled_Labels,Names_Conditions,ColorState,fs);
% SAVE ENSEMBLES DATA
% >>[Features_Ensemble,Features_Condition]=get_ensembles_features(R_Condition,Ensemble_Threshold,Ensembled_Labels,fs);
% >>save_features_ensembles(Experiment,Names_Conditions,Features_Ensemble,Features_Condition)

% ENSEMBLES FEATURES BOXPLOTS
% >>Ensembles_Features_Display

% COLOCALIZATION OF MARKED CELLS
% >>Merge_Finder_Magic

% INPSECTION to RETRIEVE ORIGINAL SIGNALS from RASTER SELECTION
% 1) Plot Raster (without sorting) from:
% >>[Rsel,IndexSorted]=Retrieve_Selected_Signal(Onsets,R_Condition,RASTER,XY_merged,XY);
% Each Rsel{n} is the selected raster
% >>Rjunto=[Rsel{1},Rsel{2},Rsel{3}]; Plot_Raster_V(Rjunto);
% Find Cell of Interest: Ci
% >>[XS,IndexSorted]=Retrieve_Selected_Signal(Onsets,R_Condition,SIGNALSclean,XY_merged,XY);
% >>figure; plot(XS{c}(Ci,:))

% LOAD RASTER FEATURES FOR A SET OF EXPERIMENTS
% >> Raster_Features_Display

%% Legacy Code:
%   Manual_Driver_Raster_Magic_Ultimate (save as private mine)