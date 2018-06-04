%% NOTE A:
% Manual Mode->Necessary to know the statistics power
% of the automatic method by dividing in -+ and --
%% RIGHT NOW


% Threshold to get NETWORK !!!!!!!!
% Add GrayChannel to modify contrast to MERGE Colocalization
% Check Driver after processing-> some +- accepted



% NO TO SO URGENT
% Re consider Accepted and Rejected ones @ automatic mode
% Get Raster Mode -> update in Manual Mode
% Driver Analysis-> Consider Derivative or Valleys
% Manual Mode for a specific raster ONLY!
% Manual mode without using pause---update workspace automatically

%% STEPS ******************************************
% ACTUAL MODE
% [RASTER_Selected_Clean,XY_selected,R_Condition,Onsets]= Select_Raster_for_NN(fs,Raster_Condition,XY,Names_Conditions,Experiment);
% R=RASTER_Selected_Clean';
% R_CONDTION1=R_Condition{1}';
% ...
% R_CONDTIONi=R_Condition{i}';
% 'Got to Dir' 
% C:\Users\Vladimir\Documents\Doctorado\Software\GetTransitum\Calcium Imaging Signal Processing\NeuralNetworks
% NeuralNetwork-> GUI mode-> Clustering Analysis
% Ensemble_Sorting
% [XY_merged,ColocateIndx]=get_merged_coordinates(Experiment,XY_selected,r);
% Plot_Ensembles_Experiment(R_Condition,EnsembleName,Ensembled_Labels,Ensemble_Threshold,UniRMutiE,ColorState,fs,[]);
% [Features_Ensemble,Features_Condition]=get_ensembles_features(R_Condition,Ensemble_Threshold,Ensembled_Labels,fs);


%% FUTURE **********************************
% Figure: reason whi mean(ROI) withput distortion
% Load Raw FLuorescenc vs F_0 distortion

% 1st Part Automatic: Raster Method
% Analyze Rejects Ones Anyway to infer Artifacts

% Set Default Settings Script / Set Default Settings (mine)



% CLUSTERING STUFF ***********************
% Threshold: prior numbercoactivyt:
% [THCOAC]=mctest(R,'modes')
% Clustering
% [Ensembles,N_Ensembles,MethodsClustering,ThEffective,
% label_frames,
% signifi_frames]=ensemble_clusterin(R,THCOAC);

% Still some Detrending issues: final points distortion !!!!!
% Initial distoriton creates valleys [SOLVED]
% (really) slow peaks [STILL]->delete usgin drive and convolution recon
% detection of modes: starter or later mode  comparisson

% Setup Intel/Info .mat File-> Default User Direcotry to save info
% Setup Script



%%CORRELATION***********************************************
simindex=corrcoef(RASTER_Selected_Clean');
cluster_index=clusterModularity(simindex,1000);
plotClusterRaster(RASTER_Selected_Clean,cluster_index,1);