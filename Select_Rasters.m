% Select_Raster Script
% Script to Run and I t Creates Variable & Variable Names
% Without User Intervention
Update_Directory;
[RASTER_Selected_Clean,XY_selected,R_Condition,Onsets]= Select_Raster_for_NN(fs,Raster_Condition,XY,Names_Conditions,Experiment);