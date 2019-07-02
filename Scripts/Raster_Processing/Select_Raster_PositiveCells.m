%% Select_Raster Script and RASTER FEATURES for POSITIVE CELLS***************
% INPUT:
%   MetaDataColocaliation.PositiveCells: Indexes of the POSITIVE SELECTED
%   fs:                 Sampling Frequency
%   R_merged:           Cell off POSITIVE SELECTED cells x Frames
%   XY_merged:          POSITIVE SELECTED Coordinates
%   Names_Conditions:   Condition Names
%   Experiment:         Experiment ID
%   CleanSignals:       POSITIVE SELECTED Denoised Signals (Transients Computing)
% Output
% @ mat File: DOESN'T UPDATE
% Features @ CSV File in ../Raster Features 
%% DO THE STUFF ***********************************************************
Update_Directory;
% Get Coordinates
XY_merged=XY_selected(MetaDataColocaliation.PositiveCells,:);
% Indexes of the Original Denoised Signals and Vectors
[CleanSignals,IndexSorted]=Retrieve_Selected_Signal(Onsets,R_Condition,ESTSIGNALS,XY_merged,XY);
% GUI-> Set OK to the Interval Times
Select_Raster_for_NN(fs,R_merged,XY_merged,Names_Conditions,Experiment,CleanSignals,Onsets);

%%  NOTES:
% Is XY.merged the same as XY(IndexSorted,:) -> HELLYES!!
% Is CleanSignals{k} the same as ESTSIGNALS{:,k}(IndexSorted,Onset:Duration) -> HELLYES!!
% Are CleanSignals{k} elicited from the Driver in R_merged{x}: -> HELLYES!


%% END OF THE WORLD