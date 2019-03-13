
%% PLOTS
[New_Index_Merged,~,RASTER_COND]=SortNeuronsCondition(R_merged);
% Plot_Raster_V(RASTER_WHOLE(New_Index,:),fs);
RASTER_COND=RASTER_COND(New_Index_Merged,:);

Plot_Raster_Ensembles(RASTER_COND,fs)
Label_Condition_Raster(Names_Conditions,R_Condition,fs);   % Labels
set(gcf,'Name',['ID: ',Experiment(2:end),' ',MetaDataColocaliation.Cells{1},'+ '],'NumberTitle','off')

[New_Index_Nomerged,~,RASTER_COND]=SortNeuronsCondition(R_nomerged);
% Plot_Raster_V(RASTER_WHOLE(New_Index,:),fs);
RASTER_COND=RASTER_COND(New_Index_Nomerged,:);

Plot_Raster_Ensembles(RASTER_COND,fs)
Label_Condition_Raster(Names_Conditions,R_Condition,fs);   % Labels
set(gcf,'Name',['ID: ',Experiment(2:end),' ',MetaDataColocaliation.Cells{1},'- '],'NumberTitle','off')
