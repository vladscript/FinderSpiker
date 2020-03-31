%% Select_Raster Script and RASTER FEATURES ********************************
% INPUT:
%   fs,
%   Raster_Condition:   Cell off ALL cells x Frames
%   XY:                 ALL Coordinates
%   Names_Conditions:   Condition Names
%   Experiment:         Experiment ID
%   ESTSIGNALS:         Denoised Signals (Transients Computing)
% Output
% @ mat File: Raster Selected, Onsets:
%   RASTER_Selected_Clean
%   XY_selected 
%   R_Condition
%   Onsets
%   New_Index_Sel
% Features @ CSV File in ../Raster Features
%% DO THE STUFF
Import_FinderSpiker;
% Raster fFeatures: denoised o sparse signal ?-> DENOISED
% Sparse Signal got a lot of peaks and false transients
Select_Raster_for_NN(fs,Raster_Condition,XY,Names_Conditions,Experiment,ESTSIGNALS);
ActualDir=pwd;
Slashes=find(ActualDir=='\');
RasterFeatsDir=[ActualDir(1:Slashes(end)),'Raster Features'];
fprintf('\n\n Move to a proper directory the CSV files of Raster features from:\n')
fprintf('%s\n\n',RasterFeatsDir);
clear;
fprintf('>>For further analysis reload the .mat file.\n')





%% END OF THE WORLD