% Script to Load & Process Calcium Fluorescence Signals (ImPatch Experiments)
% It detects Calcium Transients from VIDEOS of a single slice (EXPERIMENT)
% Sugggested Directorie's Strcuture:
% NAME_EXPERIMENT/ {VIDEOS & Coordinates from fSIENN):
% Input (same directory): 
% Coordinates:
%   ALL (XY).csv        Coordinates from 4th row x,y,r
%   or
%   ROIs (ImageJ file)
% Video Data:
%   Condition_A_00.avi
%   Condition_A_01.avi
%   ...
%   Condition_C_00.avi
%   ... 
%   Condition_Z_##.avi
% Outputs
%   ExperimentID.mat @ cd ..\Processed Data: Useful Workspace Variables
%   ExperimentID.csv @ cd ..\Features Tables:Processing Signal Features
%   ExperimentID.csv @ cd ..\Resume Tables:  Experiment Resume Features

%% Global Setup ***********************************************************
InitGlobalVars;
Load_Default_Directories;
%% Read Names, Path and Coordinates***********************************
[Names_Conditions,NumberofVideos,FN,PathName,XY,r]=Read_Videos(DefaultPath);

%% Read Experiment Videos and Sampling Frequency & Metadata
[fs,dyename,NV,NC]=ExperimentSize(NumberofVideos,FN);

%% Output Variables
InitOutVars;

%% Experiment ID
Experiment=getExpID(PathName);

%% Load Fluorescence Signals
[SIGNALS,r,H,W,RADroi]=FluorescenceSignals(NC,NumberofVideos,FN,PathName,XY,r);

%% Save(1/3) RAW Data 
save_rawdata(FolderNamePD,NumberofVideos,XY,r,PathName,RADroi,dyename,FN);

%% SETUP PROCESSING PARAMETERS ******************************
Load_Default_Values_SP;

%% AUTOMATIC DENOISE DETREND & DECONVOLUTION
[ESTSIGNALS,SNRwavelet,SNRindx,TAUSall,isSIGNALS,notSIGNALS]=AutoThreeD(NumberofVideos,L,p,taus_0,dyename,FolderNameRT);

%% SAVING(2/3) Automatic Processed Data  
save_processeddata(FolderNamePD,ESTSIGNALS,SNRwavelet,TAUSall,isSIGNALS,notSIGNALS);

%% Sort & Clean Rasters ***************************************************
[New_Index,Raster_Condition,RASTER_WHOLE_Clean,XY_clean]=ActivityMatrixSorted(XY);

%% SAVE (3/3 ) RESULTS
save_activitydata(FolderNamePD,New_Index,Raster_Condition,RASTER_WHOLE_Clean,XY_clean);

%% Visual Inpspection & Manual Processing ********************************* GREAT!
VisualInspectionInstructions;
%% END OF THE WORLD**************************************************   