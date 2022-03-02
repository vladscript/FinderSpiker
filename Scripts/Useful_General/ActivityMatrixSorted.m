function [New_Index,Raster_Condition,RASTER_WHOLE_Clean,XY_clean]=ActivityMatrixSorted(XY)
% Call global
global RASTER;
global Experiment;
global fs;
global Names_Conditions;
% make it nested function--->
% Sort by Activation in each Condition:
[New_Index,Raster_Condition,RASTER_WHOLE]=SortNeuronsCondition(RASTER,numel(XY(:,1)));
% Plot_Raster_V(RASTER_WHOLE(New_Index,:),fs);
RASTER_WHOLE_Clean=RASTER_WHOLE(New_Index,:);
XY_clean=XY(New_Index,:);
% Clean Raster and Coordinates
ActiveNeurons=find(sum(RASTER_WHOLE_Clean,2)>0);                % INDEX of Active NEURONS only
RASTER_WHOLE_Clean=RASTER_WHOLE_Clean(ActiveNeurons,:);
XY_clean=XY_clean(ActiveNeurons,:);   
%% PLOT RESULTS
Plot_Raster_Ensembles(RASTER_WHOLE_Clean,fs);                           % Clean Whole Raster
CurrentFig=gcf;
CurrentFig.Name = ['ID: ',Experiment,' automatically-processed'];
CurrentFig.NumberTitle='off';
Label_Condition_Raster(Names_Conditions,Raster_Condition,fs);   % Labels
