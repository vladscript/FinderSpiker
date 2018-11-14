%% NOTES##################################################################
% Manual Mode->Necessary to know the statistics power
% of the automatic method by dividing in -+ and -- (false+ & false-)
% Old Version Manual Mode:
%   Manual_Driver_Raster_Magic.m
%% FIXED  RECENTLY @ GIT
% Automatize Clustering : OK
% Bug at Experiment Name

%% Bugs & New Functions NOW

% Add button to save Zoom image (MERGED MAGIC)
% Save Points SELECTED-> add to file .mat
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
% Kalman Filtering at SNR and lambdas pdf's: for optimal threshold
% Automatize MERGE SELECTOR

% Setup Intel/Info .mat File-> Default User Direcotry to save info
% Setup Script

%% STEPS ******************************************
%PROCESSING
% >>FInder_Spiker_Calcium
% >>Detected_Visual_Inspection
% >>Undetected_Visual_Inspection
% >>Save_and_Plot

% RASTER SELECTION
% ACTUAL MODE: @ Original Coordiantes Order
% >>[RASTER_Selected_Clean,XY_selected,R_Condition,Onsets]= Select_Raster_for_NN(fs,Raster_Condition,XY,Names_Conditions,Experiment);


% >>R=RASTER_Selected_Clean';
% >>R_CONDTION1=R_Condition{1}';
% ...
% >>R_CONDTIONi=R_Condition{i}';

% LOAD RASTER FEATURES
% >> Raster_Features_Display

% CLUSTERING
% For each R_CONDITIONi:
% AUTOMATIC:
% >>R_CONDITIONi_Analysis=get_ensembles(R_CONDITIONi); 
% MANUAL
% Set CAG_TH,  Nensembles & SimMethod='hamming'
% >>Ri_Analysis = raster_cluster(R_CONDITIONi,CAG_TH,Nensembles,SimMethod); 
% DISPLAY AND SAVE RESULTS
% >> Ensemble_Sorting

% >>Plot_Ensembles_Experiment(R_Condition,EnsembleName,Ensembled_Labels,Ensemble_Threshold,UniRMutiE,ColorState,fs,[]);
% >>Plot_Hebbian_Paths(R_Condition,Ensemble_Threshold,Ensembled_Labels,Names_Conditions,ColorState,fs);
% >>[Features_Ensemble,Features_Condition]=get_ensembles_features(R_Condition,Ensemble_Threshold,Ensembled_Labels,fs);
% >>save_features_ensembles(Experiment,Names_Conditions,Features_Ensemble,Features_Condition)

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


%% Legacy Code:
%   Manual_Driver_Raster_Magic_Ultimate (save as private mine)