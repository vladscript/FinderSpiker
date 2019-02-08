%% GET RASTER
[New_Index,Raster_Condition,RASTER_WHOLE]=SortNeuronsCondition(RASTER);
RASTER_WHOLE_Clean=RASTER_WHOLE(New_Index,:);
ActiveNeurons=find(sum(RASTER_WHOLE_Clean,2)>0);                % INDEX of Active NEURONS only
RASTER_WHOLE_Clean=RASTER_WHOLE_Clean(ActiveNeurons,:);
%% PLOT RESULTS
Plot_Raster_V(RASTER_WHOLE_Clean,fs);                           % Clean Whole Raster
set(gcf,'Name',['ID: ',Experiment(2:end),' processed'],'NumberTitle','off')
Label_Condition_Raster(Names_Conditions,Raster_Condition,fs);   % Labels