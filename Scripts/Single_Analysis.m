% Script to Process Calcium Fluorescence Signals (CSV records)
% It detects Calcium Transients from CSV data of a single slice (EXPERIMENT)
% Input (same directory): 
% Coordinates:
%   ALL (XY).csv        Coordinates from 4th row x,y,r
%   or
%   ROIs (ImageJ file)
% Fluorescence Data:
%   CONDITION_A.csv
%   ...
%   Condition_K.csv
%   ... 
%   Condition_Z.csv
% Outputs
%   ExperimentID.mat @ cd ..\Processed Data: Useful Workspace Variables
%   ExperimentID.csv @ cd ..\Features Tables:Processing Signal Features
%   ExperimentID.csv @ cd ..\Resume Tables:  Experiment Resume Features
%% Global Setup ***********************************************************
InitGlobalVars;
Load_Default_Directories;

%% Read Name, Path and Coordinates (if are there)**************************
[Names_Conditions,NumberofVideos,FN,PathName,XY,r]=Read_CSV(DefaultPath);
NC=1;
NV= str2double( NumberofVideos{1} );
% Read Fluorophore DYe
dyename = inputdlg('Fluorophore : ','DYE', [1 75]);
%% Output Variables
InitOutVars; 
%% Experiment ID
Experiment=getExpID(FN{1});
%% Load Fluorescence Signals
[SIGNALS,r,H,W,RADroi,fs]=FluorescenceSignalsCSV(NV,FN,PathName,r);
stopproc=false;
if isempty(r)
    XY=zeros(size(SIGNALS{1},1),2);
    r=zeros(size(SIGNALS{1},1),1);
    fprintf('\n>[WARNING]: no coordinates provided\n')
else
    if numel(r)~=size(SIGNALS{1},1)
        warning('Missmatch of number of ROIs and Fluorescence signals:')
        fprintf('>#ROIs:%i and #Cells:%i\n>>Analysis Stopped\n',numel(r),size(SIGNALS{1},1))
        stopproc=true;
    end
end
if ~stopproc
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
end
    %% END OF THE WORLD**************************************************   