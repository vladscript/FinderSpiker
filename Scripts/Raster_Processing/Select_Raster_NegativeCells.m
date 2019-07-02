% Select_Raster Script for POSIVE CELLS
% 
Update_Directory;

% Get Coordinates
XY_nomerged=XY_selected(MetaDataColocaliation.NegativeCells,:);

% Indexes of the Original Denoised Signals and Vectors
[CleanSignals,IndexSorted]=Retrieve_Selected_Signal(Onsets,R_Condition,ESTSIGNALS,XY_nomerged,XY);

Select_Raster_for_NN(fs,R_nomerged,XY_nomerged,Names_Conditions,Experiment,CleanSignals,Onsets);


% CHECKED: OK

%% END OF THE WORLD