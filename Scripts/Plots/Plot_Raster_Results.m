%% PLOT RESULTS
Plot_Raster_V(RASTER_WHOLE_Clean,fs);                           % Clean Whole Raster
drawnow;
hfig=gcf;
set(hfig,'Name',['ID: ',Experiment(2:end)],'NumberTitle','off')
Label_Condition_Raster(Names_Conditions,Raster_Condition,fs);   % Labels